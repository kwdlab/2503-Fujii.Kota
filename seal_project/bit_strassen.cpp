#include "seal/seal.h"
#include <iostream>
#include <vector>
#include <iomanip> // chronoライブラリのインクルード


using namespace std;
using namespace seal;

vector<Ciphertext> encrypted_matrix_multiplication_strassen(
    const vector<Ciphertext> &strassen1, const vector<Ciphertext> &strassen2,
    Evaluator &evaluator, RelinKeys &relin_keys, size_t n)
{
    // ベースケース（行列サイズが1×1の場合）
    if (n == 1)
    {
        vector<Ciphertext> result(1);
        evaluator.multiply(strassen1[0], strassen2[0], result[0]);
        evaluator.relinearize_inplace(result[0], relin_keys);
        return result;
    }

    // 行列サイズを2で割った半分のサイズ
    size_t new_size = n / 2;

    // 行列の分割
    vector<Ciphertext> A11(new_size * new_size), A12(new_size * new_size),
        A21(new_size * new_size), A22(new_size * new_size);
    vector<Ciphertext> B11(new_size * new_size), B12(new_size * new_size),
        B21(new_size * new_size), B22(new_size * new_size);

    for (size_t i = 0; i < new_size; ++i)
    {
        for (size_t j = 0; j < new_size; ++j)
        {
            A11[i * new_size + j] = strassen1[i * n + j];
            A12[i * new_size + j] = strassen1[i * n + (j + new_size)];
            A21[i * new_size + j] = strassen1[(i + new_size) * n + j];
            A22[i * new_size + j] = strassen1[(i + new_size) * n + (j + new_size)];

            B11[i * new_size + j] = strassen2[i * n + j];
            B12[i * new_size + j] = strassen2[i * n + (j + new_size)];
            B21[i * new_size + j] = strassen2[(i + new_size) * n + j];
            B22[i * new_size + j] = strassen2[(i + new_size) * n + (j + new_size)];
        }
    }

    // Strassenアルゴリズムの中間計算
    vector<Ciphertext> P1(new_size * new_size), P2(new_size * new_size),
        P3(new_size * new_size), P4(new_size * new_size),
        P5(new_size * new_size), P6(new_size * new_size),
        P7(new_size * new_size);

    vector<Ciphertext> temp1(new_size * new_size), temp2(new_size * new_size);

    // P1 = (A11 + A22) * (B11 + B22)
    for (size_t i = 0; i < new_size * new_size; ++i)
    {
        evaluator.add(A11[i], A22[i], temp1[i]);
        evaluator.add(B11[i], B22[i], temp2[i]);
    }
    P1 = encrypted_matrix_multiplication_strassen(temp1, temp2, evaluator, relin_keys, new_size);

    // P2 = (A21 + A22) * B11
    for (size_t i = 0; i < new_size * new_size; ++i)
    {
        evaluator.add(A21[i], A22[i], temp1[i]);
    }
    P2 = encrypted_matrix_multiplication_strassen(temp1, B11, evaluator, relin_keys, new_size);

    // P3 = A11 * (B12 - B22)
    for (size_t i = 0; i < new_size * new_size; ++i)
    {
        evaluator.sub(B12[i], B22[i], temp1[i]);
    }
    P3 = encrypted_matrix_multiplication_strassen(A11, temp1, evaluator, relin_keys, new_size);

    // P4 = A22 * (B21 - B11)
    for (size_t i = 0; i < new_size * new_size; ++i)
    {
        evaluator.sub(B21[i], B11[i], temp1[i]);
    }
    P4 = encrypted_matrix_multiplication_strassen(A22, temp1, evaluator, relin_keys, new_size);

    // P5 = (A11 + A12) * B22
    for (size_t i = 0; i < new_size * new_size; ++i)
    {
        evaluator.add(A11[i], A12[i], temp1[i]);
    }
    P5 = encrypted_matrix_multiplication_strassen(temp1, B22, evaluator, relin_keys, new_size);

    // P6 = (A21 - A11) * (B11 + B12)
    for (size_t i = 0; i < new_size * new_size; ++i)
    {
        evaluator.sub(A21[i], A11[i], temp1[i]);
        evaluator.add(B11[i], B12[i], temp2[i]);
    }
    P6 = encrypted_matrix_multiplication_strassen(temp1, temp2, evaluator, relin_keys, new_size);

    // P7 = (A12 - A22) * (B21 + B22)
    for (size_t i = 0; i < new_size * new_size; ++i)
    {
        evaluator.sub(A12[i], A22[i], temp1[i]);
        evaluator.add(B21[i], B22[i], temp2[i]);
    }
    P7 = encrypted_matrix_multiplication_strassen(temp1, temp2, evaluator, relin_keys, new_size);

    // 部分行列の統合
    vector<Ciphertext> C11(new_size * new_size), C12(new_size * new_size),
        C21(new_size * new_size), C22(new_size * new_size);

    for (size_t i = 0; i < new_size * new_size; ++i)
    {
        // C11 = P1 + P4 - P5 + P7
        evaluator.add(P1[i], P4[i], C11[i]);
        evaluator.sub_inplace(C11[i], P5[i]);
        evaluator.add_inplace(C11[i], P7[i]);

        // C12 = P3 + P5
        evaluator.add(P3[i], P5[i], C12[i]);

        // C21 = P2 + P4
        evaluator.add(P2[i], P4[i], C21[i]);

        // C22 = P1 + P3 - P2 + P6
        evaluator.add(P1[i], P3[i], C22[i]);
        evaluator.sub_inplace(C22[i], P2[i]);
        evaluator.add_inplace(C22[i], P6[i]);
    }

    // 結果の結合
    vector<Ciphertext> result(n * n);
    for (size_t i = 0; i < new_size; ++i)
    {
        for (size_t j = 0; j < new_size; ++j)
        {
            result[i * n + j] = C11[i * new_size + j];
            result[i * n + (j + new_size)] = C12[i * new_size + j];
            result[(i + new_size) * n + j] = C21[i * new_size + j];
            result[(i + new_size) * n + (j + new_size)] = C22[i * new_size + j];
        }
    }

    return result;
}


// 結果を復号化して表示する関数
void display_matrix_result_batch(const vector<Ciphertext> &encrypted_matrix,
                                 Decryptor &decryptor, BatchEncoder &batch_encoder, size_t n) {
    vector<uint64_t> decoded_result(batch_encoder.slot_count());

    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            Plaintext plain_result;
            decryptor.decrypt(encrypted_matrix[i * n + j], plain_result);
            batch_encoder.decode(plain_result, decoded_result);
            cout << decoded_result[0] << " "; // スロット0の値を表示
        }
        cout << endl;
    }
}

// 平文行列の初期化を動的に行う
vector<uint64_t> generate_plain_matrix(size_t n, bool random = false) {
    vector<uint64_t> matrix(n * n);
    for (size_t i = 0; i < n; ++i) {
        for (size_t j = 0; j < n; ++j) {
            matrix[i * n + j] = random ? rand() % 10 + 1 : (i * n + j + 1); // 1から10のランダム値またはデフォルト値
        }
    }
    return matrix;
}

int main() {
    // SEALパラメータの設定
    EncryptionParameters parms(scheme_type::bfv);
    size_t poly_modulus_degree = 16384;
    parms.set_poly_modulus_degree(poly_modulus_degree);
    parms.set_coeff_modulus(CoeffModulus::BFVDefault(poly_modulus_degree));
    parms.set_plain_modulus(PlainModulus::Batching(poly_modulus_degree, 20)); // ビット設定

    SEALContext context(parms);

    // キーの生成
    KeyGenerator keygen(context);
    SecretKey secret_key = keygen.secret_key();
    PublicKey public_key;
    keygen.create_public_key(public_key);
    RelinKeys relin_keys;
    keygen.create_relin_keys(relin_keys);

    Encryptor encryptor(context, public_key);
    Decryptor decryptor(context, secret_key);
    Evaluator evaluator(context);
    BatchEncoder batch_encoder(context);

    // 行列のサイズ（n×n）: 動的に設定可能
    size_t n = 4; // 行列サイズ
    
    //行列文字（trueならランダム、falseならシーケンス（１、２、３...））
    vector<uint64_t> strassen1_plain = generate_plain_matrix(16, false);
    vector<uint64_t> strassen2_plain = generate_plain_matrix(16, false);

    
    // *** 実行開始からストラッセン法終了までの計測開始 ***
    auto start_time = chrono::high_resolution_clock::now();
    
    // 平文行列を暗号化
    vector<Ciphertext> strassen1_encrypted(n * n);
    vector<Ciphertext> strassen2_encrypted(n * n);

    for (size_t i = 0; i < n * n; ++i) {
        Plaintext plain1, plain2;
        batch_encoder.encode(vector<uint64_t>(batch_encoder.slot_count(), strassen1_plain[i]), plain1);
        batch_encoder.encode(vector<uint64_t>(batch_encoder.slot_count(), strassen2_plain[i]), plain2);

        encryptor.encrypt(plain1, strassen1_encrypted[i]);
        encryptor.encrypt(plain2, strassen2_encrypted[i]);
    }

    /// 暗号化された行列の掛け算
    vector<Ciphertext> encrypted_result = encrypted_matrix_multiplication_strassen(
        strassen1_encrypted, strassen2_encrypted, evaluator, relin_keys, n);
    
    // *** 計測終了 ***
    auto end_time = chrono::high_resolution_clock::now();
    
    // 実行時間の表示
    cout << fixed << setprecision(3)
         << "Execution Time (Start to End of Strassen Calculation): "
         << chrono::duration<double>(end_time - start_time).count() << " s" << endl;

    return 0;
}
